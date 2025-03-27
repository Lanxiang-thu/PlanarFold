#!/bin/bash

# 提示用户输入 Conda 环境名称
read -p "请输入您的 Conda 环境名称（默认为 my_python3.9）: " new_shebang


# 定义旧 Shebang 路径（需与你的原路径一致）
old_shebang="#!/home/xiangl/miniconda3/envs/my_python3.9/bin/python"


# 遍历所有文件并替换 Shebang
find "../tools/" -type f \( -name "*.py" -o -name "*.sh" \) -print0 | while IFS= read -r -d $'\0' file; do
  # 检查第一行是否为旧 Shebang
  if head -n 1 "$file" | grep -qF "$old_shebang"; then
    echo "修改文件: $file"
    # 根据系统类型处理 sed 命令（兼容 macOS 和 Linux）
    if sed --version 2>/dev/null | grep -q "GNU"; then
      sed -i "1s|${old_shebang}|#!${new_shebang}|" "$file"
    else
      sed -i '' "1s|${old_shebang}|#!${new_shebang}|" "$file"
    fi
  fi
done

echo "Replacement finished: ${new_shebang}"
